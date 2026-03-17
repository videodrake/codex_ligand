#!/usr/bin/env python3
"""Phase 1 Input Preparation: Full kinase domain + extended beta-meander.

Task Group 1.0 implementation for Phase 1 PPI-First Interface Mapping (v2).

This script prepares structurally adequate receptor and partner inputs that
eliminate known fragment-docking artifacts:

  1. Full kinase domain receptor (~280-309 residues, N-lobe + C-lobe)
     for each of the three receptor states (3GT8_raw, cl38_48, cl85_100)
  2. Extended beta-meander partner (~955-1006, ~52 residues)
     extracted from the TH1 domain structure
  3. Docking-ready 2-chain PDBs (receptor=A, partner=B)
  4. Metadata, validation report, and pilot data registration

Numbering system: All residue numbers are PDB numbering (3GT8-consistent).
3GT8 PDB numbering differs from UniProt (P00533) by approximately +24.

Usage:
    python -m egfr_pipeline.phase1.prepare_inputs [--output_dir INPUT_DIR]
"""

import argparse
import csv
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE1_RUNTIME_INPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi" / "runtime_inputs"


def _display_path(path: Path) -> str:
    """Return a project-relative path when possible, else an absolute path."""
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Receptor states and their source PDB files
RECEPTOR_STATES = {
    "3GT8_raw": {
        "pdb": "input/receptors/3GT8_raw.pdb",
        "source_chain": "A",       # Crystal structure, use chain A monomer
        "description": "Crystal structure (PDB 3GT8), chain A monomer",
    },
    "EGFR_160-185": {
        "pdb": "input/receptors/EGFR_160-185.pdb",
        "source_chain": "X",       # MD cluster, single chain X
        "description": "MD cluster representative (38-48 ns), single chain",
    },
    "EGFR_170-200": {
        "pdb": "input/receptors/EGFR_170-200.pdb",
        "source_chain": "X",       # MD cluster, single chain X
        "description": "MD cluster representative (85-100 ns), single chain",
    },
}

# Common residue range for kinase domain (PDB numbering)
# 3GT8_raw chain A: 699-1007 (309 residues)
# MD clusters chain X: 634-1014 (381 residues, includes flanking)
# Overlap for comparability: 699-1007
KINASE_DOMAIN_START = 699
KINASE_DOMAIN_END = 1007

# Approximate N-lobe / C-lobe boundary (PDB numbering)
# The hinge region is approximately at residue 838 (PDB numbering)
NLOBE_CLOBE_BOUNDARY = 838

# Extended beta-meander definition
BETA_MEANDER_SOURCE = "input/PPI/TH1 domain.pdb"
BETA_MEANDER_SOURCE_CHAIN = "A"
BETA_MEANDER_START = 955   # Extended from original 960/962
BETA_MEANDER_END = 1006

# Sheet definitions (MYO1D numbering = PDB numbering for this structure)
SHEET_DEFINITIONS = {
    "sheet_8": {"residues": [961, 962, 963, 964], "role": "active_face_primary"},
    "sheet_9": {"residues": [968, 969, 970, 971, 972], "role": "active_face_primary"},
    "sheet_10": {"residues": [977, 978, 979, 980], "role": "neutral"},
    "sheet_11": {"residues": [985, 986, 987, 988], "role": "neutral"},
    "sheet_12": {"residues": [993, 994, 995, 996, 997], "role": "structural_support"},
}

# Membrane-proximal residues to exclude (PDB numbering, monomer chain A only)
MEMBRANE_PROXIMAL = "709-720,724-731,736-739,747,783-785,799-805,871-873,917-921"

# Pilot data reference (legacy C-lobe fragment system)
PILOT_DATA = {
    "construct_type": "legacy_clobe_fragment",
    "receptor_residues": "960-1006",  # ~45 residues C-lobe fragment
    "partner_residues": "960-1006",   # truncated beta-meander
    "partner_start": 960,
    "known_artifacts": [
        "N-lobe absence: C-lobe surface landscape differs from full kinase domain",
        "VAL962 N-terminal artifact: first residue has artificial charge/freedom",
        "No orientation filtering applied",
    ],
}


# ---------------------------------------------------------------------------
# PDB manipulation helpers
# ---------------------------------------------------------------------------

def parse_pdb_atoms(pdb_path: Path) -> List[str]:
    """Read PDB and return ATOM lines (protein only)."""
    lines = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM"):
                lines.append(line)
    return lines


def get_chain(line: str) -> str:
    """Extract chain ID from PDB ATOM line."""
    return line[21]


def get_resnum(line: str) -> int:
    """Extract residue number from PDB ATOM line."""
    return int(line[22:26].strip())


def get_resname(line: str) -> str:
    """Extract residue name from PDB ATOM line."""
    return line[17:20].strip()


def get_atom_name(line: str) -> str:
    """Extract atom name from PDB ATOM line."""
    return line[12:16].strip()


def set_chain(line: str, chain: str) -> str:
    """Set chain ID in PDB ATOM line."""
    return line[:21] + chain + line[22:]


def set_atom_serial(line: str, serial: int) -> str:
    """Set atom serial number in PDB ATOM line."""
    return line[:6] + f"{serial:5d}" + line[11:]


# Required backbone atoms for PyRosetta pose loading
_BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def check_backbone_completeness(
    atom_lines: List[str],
) -> Tuple[List[int], Dict[int, Set[str]]]:
    """Check backbone atom completeness for each residue.

    Returns (incomplete_resnums, missing_map) where missing_map maps
    resnum -> set of missing backbone atom names.
    """
    residue_atoms: Dict[int, Set[str]] = {}
    for line in atom_lines:
        resnum = get_resnum(line)
        residue_atoms.setdefault(resnum, set()).add(get_atom_name(line))

    incomplete: List[int] = []
    missing_map: Dict[int, Set[str]] = {}
    for resnum in sorted(residue_atoms):
        missing = _BACKBONE_ATOMS - residue_atoms[resnum]
        if missing:
            incomplete.append(resnum)
            missing_map[resnum] = missing

    return incomplete, missing_map


def strip_incomplete_residues(
    atom_lines: List[str],
    resnums: List[int],
    incomplete_resnums: List[int],
) -> Tuple[List[str], List[int], List[int]]:
    """Remove residues with incomplete backbone atoms.

    Returns (cleaned_lines, cleaned_resnums, removed_resnums).
    """
    remove_set = set(incomplete_resnums)
    cleaned_lines = [l for l in atom_lines if get_resnum(l) not in remove_set]
    cleaned_resnums = [r for r in resnums if r not in remove_set]
    return cleaned_lines, cleaned_resnums, sorted(remove_set)


def extract_chain_residue_range(
    pdb_path: Path,
    source_chain: str,
    res_start: int,
    res_end: int,
    target_chain: str = "A",
) -> Tuple[List[str], List[int]]:
    """Extract atoms from a specific chain and residue range.

    Returns (atom_lines, unique_residue_numbers).
    """
    all_atoms = parse_pdb_atoms(pdb_path)
    selected = []
    resnums = set()

    for line in all_atoms:
        chain = get_chain(line)
        resnum = get_resnum(line)
        if chain == source_chain and res_start <= resnum <= res_end:
            if chain != target_chain:
                line = set_chain(line, target_chain)
            selected.append(line)
            resnums.add(resnum)

    return selected, sorted(resnums)


def write_pdb(lines: List[str], output_path: Path, renumber_atoms: bool = True) -> None:
    """Write PDB file with proper TER/END records and atom renumbering."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        prev_chain = None
        serial = 0
        for line in lines:
            serial += 1
            if renumber_atoms:
                line = set_atom_serial(line, serial)
            curr_chain = get_chain(line)
            if prev_chain is not None and curr_chain != prev_chain:
                f.write("TER\n")
            f.write(line)
            prev_chain = curr_chain
        f.write("TER\nEND\n")


def get_residue_sequence(pdb_path: Path, chain: str) -> List[Tuple[int, str]]:
    """Get ordered list of (resnum, resname) for a chain."""
    atoms = parse_pdb_atoms(pdb_path)
    seen = {}
    for line in atoms:
        if get_chain(line) == chain:
            rn = get_resnum(line)
            if rn not in seen:
                seen[rn] = get_resname(line)
    return [(rn, seen[rn]) for rn in sorted(seen)]


def parse_ranges(range_str: str) -> Set[int]:
    """Parse '699-711,748-770' into set of ints."""
    result = set()
    for part in range_str.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            result.update(range(int(a), int(b) + 1))
        else:
            result.add(int(part))
    return result


def format_ranges(nums) -> str:
    """Format a set/list of ints into '699-701,748-749'."""
    if not nums:
        return ""
    sorted_nums = sorted(nums)
    ranges = []
    start = prev = sorted_nums[0]
    for n in sorted_nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append(f"{start}-{prev}" if start != prev else str(start))
            start = prev = n
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(ranges)


# ---------------------------------------------------------------------------
# Task 1.0.1: Full kinase domain receptor preparation
# ---------------------------------------------------------------------------

def prepare_receptor(
    state_name: str,
    state_info: dict,
    output_dir: Path,
) -> dict:
    """Prepare a full kinase domain receptor PDB.

    For 3GT8_raw: extract chain A monomer (699-1007).
    For MD clusters: extract chain X residues 699-1007, rename to chain A.

    Returns metadata dict.
    """
    pdb_path = PROJECT_ROOT / state_info["pdb"]
    source_chain = state_info["source_chain"]

    print(f"\n{'='*60}")
    print(f"Preparing receptor: {state_name}")
    print(f"  Source: {pdb_path}")
    print(f"  Source chain: {source_chain}")
    print(f"  Target range: {KINASE_DOMAIN_START}-{KINASE_DOMAIN_END} (chain A)")

    # Extract kinase domain
    lines, resnums = extract_chain_residue_range(
        pdb_path, source_chain,
        KINASE_DOMAIN_START, KINASE_DOMAIN_END,
        target_chain="A",
    )

    if not resnums:
        raise ValueError(
            f"No residues found in {pdb_path} chain {source_chain} "
            f"range {KINASE_DOMAIN_START}-{KINASE_DOMAIN_END}"
        )

    # --- Backbone completeness check ---
    incomplete, missing_bb = check_backbone_completeness(lines)
    stripped_resnums: List[int] = []
    if incomplete:
        lines, resnums, stripped_resnums = strip_incomplete_residues(
            lines, resnums, incomplete,
        )
        for rn in stripped_resnums:
            print(f"  STRIPPED residue {rn}: missing backbone {missing_bb[rn]}")
        if not resnums:
            raise ValueError(
                f"All residues in {pdb_path} had incomplete backbone — "
                f"source PDB may be corrupt"
            )

    n_residues = len(resnums)
    n_atoms = len(lines)
    first_res = resnums[0]
    last_res = resnums[-1]

    # Check for gaps
    expected = set(range(first_res, last_res + 1))
    actual = set(resnums)
    missing = sorted(expected - actual)

    # N-lobe / C-lobe counts
    n_nlobe = sum(1 for r in resnums if r < NLOBE_CLOBE_BOUNDARY)
    n_clobe = sum(1 for r in resnums if r >= NLOBE_CLOBE_BOUNDARY)

    # Write receptor PDB
    output_pdb = output_dir / f"receptor_{state_name}.pdb"
    write_pdb(lines, output_pdb)

    print(f"  Output: {output_pdb}")
    print(f"  Residues: {first_res}-{last_res} ({n_residues} residues, {n_atoms} atoms)")
    if stripped_resnums:
        print(f"  Backbone-incomplete residues removed: {format_ranges(stripped_resnums)}")
    print(f"  N-lobe: {n_nlobe} residues (< {NLOBE_CLOBE_BOUNDARY})")
    print(f"  C-lobe: {n_clobe} residues (>= {NLOBE_CLOBE_BOUNDARY})")
    if missing:
        print(f"  WARNING: {len(missing)} missing residues: {format_ranges(missing)}")
    else:
        print(f"  No gaps detected")

    # Get sequence for validation
    seq = get_residue_sequence(output_pdb, "A")

    metadata = {
        "state_name": state_name,
        "construct_type": "full_kinase_domain",
        "source_pdb": str(state_info["pdb"]),
        "source_chain": source_chain,
        "description": state_info["description"],
        "output_pdb": _display_path(output_pdb),
        "chain_id": "A",
        "residue_start": first_res,
        "residue_end": last_res,
        "n_residues": n_residues,
        "n_atoms": n_atoms,
        "n_nlobe_residues": n_nlobe,
        "n_clobe_residues": n_clobe,
        "nlobe_clobe_boundary": NLOBE_CLOBE_BOUNDARY,
        "missing_residues": format_ranges(missing) if missing else "",
        "numbering_system": "PDB (3GT8-consistent)",
        "membrane_proximal_excluded": MEMBRANE_PROXIMAL,
    }

    return metadata


# ---------------------------------------------------------------------------
# Task 1.0.3: Extended beta-meander partner preparation
# ---------------------------------------------------------------------------

def prepare_extended_beta_meander(output_dir: Path) -> dict:
    """Prepare extended beta-meander partner (~955-1006).

    Extracts from TH1 domain structure, which contains the full 801-1006 range.

    Returns metadata dict.
    """
    th1_path = PROJECT_ROOT / BETA_MEANDER_SOURCE

    print(f"\n{'='*60}")
    print(f"Preparing extended beta-meander partner")
    print(f"  Source: {th1_path}")
    print(f"  Target range: {BETA_MEANDER_START}-{BETA_MEANDER_END}")

    # Extract extended beta-meander
    lines, resnums = extract_chain_residue_range(
        th1_path, BETA_MEANDER_SOURCE_CHAIN,
        BETA_MEANDER_START, BETA_MEANDER_END,
        target_chain="A",  # standalone file uses chain A
    )

    if not resnums:
        raise ValueError(
            f"No residues found in {th1_path} chain {BETA_MEANDER_SOURCE_CHAIN} "
            f"range {BETA_MEANDER_START}-{BETA_MEANDER_END}"
        )

    # --- Backbone completeness check ---
    incomplete, missing_bb = check_backbone_completeness(lines)
    stripped_resnums: List[int] = []
    if incomplete:
        lines, resnums, stripped_resnums = strip_incomplete_residues(
            lines, resnums, incomplete,
        )
        for rn in stripped_resnums:
            print(f"  STRIPPED residue {rn}: missing backbone {missing_bb[rn]}")
        if not resnums:
            raise ValueError(
                f"All residues in {th1_path} had incomplete backbone — "
                f"source PDB may be corrupt"
            )

    n_residues = len(resnums)
    n_atoms = len(lines)
    first_res = resnums[0]
    last_res = resnums[-1]

    # Validate first residue is NOT VAL962 (the old truncation artifact)
    seq = []
    seen = set()
    for line in lines:
        rn = get_resnum(line)
        if rn not in seen:
            seen.add(rn)
            seq.append((rn, get_resname(line)))

    first_resname = seq[0][1] if seq else "?"
    val962_is_first = (first_res == 962)

    # Check gaps
    expected = set(range(first_res, last_res + 1))
    actual = set(resnums)
    missing = sorted(expected - actual)

    # Annotate sheets
    sheet_annotations = {}
    for sheet_name, info in SHEET_DEFINITIONS.items():
        present = [r for r in info["residues"] if r in actual]
        total = len(info["residues"])
        sheet_annotations[sheet_name] = {
            "residues": info["residues"],
            "role": info["role"],
            "present_count": len(present),
            "total_count": total,
            "complete": len(present) == total,
        }

    # Write standalone partner PDB
    output_pdb = output_dir / "partner_extended_beta_meander.pdb"
    write_pdb(lines, output_pdb)

    print(f"  Output: {output_pdb}")
    print(f"  Residues: {first_res}-{last_res} ({n_residues} residues, {n_atoms} atoms)")
    print(f"  First residue: {first_resname}{first_res}")
    print(f"  VAL962 is first residue: {val962_is_first}")
    if val962_is_first:
        print(f"  WARNING: VAL962 is still the first residue! Extension may have failed.")
    else:
        print(f"  OK: VAL962 is NOT the first residue (artifact eliminated)")
    if missing:
        print(f"  WARNING: {len(missing)} missing residues: {format_ranges(missing)}")

    # Print sheet completeness
    print(f"\n  Sheet annotations:")
    for sname, sinfo in sheet_annotations.items():
        status = "COMPLETE" if sinfo["complete"] else f"INCOMPLETE ({sinfo['present_count']}/{sinfo['total_count']})"
        print(f"    {sname}: {sinfo['residues']} [{sinfo['role']}] — {status}")

    metadata = {
        "construct_type": "extended_beta_meander",
        "source_pdb": BETA_MEANDER_SOURCE,
        "source_chain": BETA_MEANDER_SOURCE_CHAIN,
        "output_pdb": _display_path(output_pdb),
        "chain_id": "A",
        "residue_start": first_res,
        "residue_end": last_res,
        "n_residues": n_residues,
        "n_atoms": n_atoms,
        "first_residue": f"{first_resname}{first_res}",
        "val962_is_first": val962_is_first,
        "missing_residues": format_ranges(missing) if missing else "",
        "sheet_annotations": sheet_annotations,
        "numbering_system": "MYO1D PDB numbering",
        "active_face_residues": SHEET_DEFINITIONS["sheet_8"]["residues"] + SHEET_DEFINITIONS["sheet_9"]["residues"],
        "structural_support_residues": SHEET_DEFINITIONS["sheet_12"]["residues"],
    }

    return metadata


# ---------------------------------------------------------------------------
# Task 1.0.1 + 1.0.3 combined: Docking-ready 2-chain PDBs
# ---------------------------------------------------------------------------

def prepare_docking_pair(
    receptor_pdb: Path,
    partner_pdb: Path,
    output_pdb: Path,
    state_name: str,
) -> dict:
    """Create a docking-ready PDB with receptor=chain A, partner=chain B.

    This is the Phase 1 approach: single kinase domain monomer + partner.
    No dimer merging (unlike the legacy C-lobe fragment approach).
    """
    print(f"\n  Creating docking pair: {state_name}")

    # Read receptor (already chain A)
    receptor_lines = parse_pdb_atoms(receptor_pdb)
    r_resnums = sorted(set(get_resnum(l) for l in receptor_lines))

    # Read partner and set to chain B
    partner_lines = parse_pdb_atoms(partner_pdb)
    partner_b_lines = [set_chain(l, "B") for l in partner_lines]
    p_resnums = sorted(set(get_resnum(l) for l in partner_lines))

    # --- Final backbone completeness check on combined PDB ---
    combined = receptor_lines + partner_b_lines
    incomplete, missing_bb = check_backbone_completeness(combined)
    stripped: List[int] = []
    if incomplete:
        combined, _, stripped = strip_incomplete_residues(
            combined, r_resnums + p_resnums, incomplete,
        )
        for rn in stripped:
            print(f"    STRIPPED residue {rn}: missing backbone {missing_bb[rn]}")
        # Rebuild per-chain residue lists from cleaned lines
        r_resnums = sorted(set(
            get_resnum(l) for l in combined if get_chain(l) == "A"
        ))
        p_resnums = sorted(set(
            get_resnum(l) for l in combined if get_chain(l) == "B"
        ))

    write_pdb(combined, output_pdb)

    # Compute excluded residues
    excl_set = parse_ranges(MEMBRANE_PROXIMAL)
    # Only keep excluded residues that are actually in this receptor
    excl_actual = excl_set & set(r_resnums)
    excl_str = format_ranges(excl_actual)

    print(f"    Receptor: chain A, {r_resnums[0]}-{r_resnums[-1]} ({len(r_resnums)} res)")
    print(f"    Partner:  chain B, {p_resnums[0]}-{p_resnums[-1]} ({len(p_resnums)} res)")
    if stripped:
        print(f"    Backbone-incomplete residues removed: {format_ranges(stripped)}")
    print(f"    Excluded: {len(excl_actual)} membrane-proximal residues")
    print(f"    Output:   {output_pdb}")

    return {
        "state_name": state_name,
        "output_pdb": _display_path(output_pdb),
        "receptor_chain": "A",
        "receptor_range": f"{r_resnums[0]}-{r_resnums[-1]}",
        "receptor_n_residues": len(r_resnums),
        "partner_chain": "B",
        "partner_range": f"{p_resnums[0]}-{p_resnums[-1]}",
        "partner_n_residues": len(p_resnums),
        "excluded_residues_A": excl_str,
        "n_excluded": len(excl_actual),
        "total_atoms": len(combined),
        "stripped_residues": format_ranges(stripped) if stripped else "",
    }


# ---------------------------------------------------------------------------
# Task 1.0.2: Receptor validation
# ---------------------------------------------------------------------------

def validate_receptors(receptor_metas: List[dict]) -> List[str]:
    """Cross-validate all receptor states for comparability.

    Returns list of validation messages.
    """
    messages = []
    messages.append("=" * 60)
    messages.append("RECEPTOR CROSS-VALIDATION")
    messages.append("=" * 60)

    # Check chain ID consistency
    chains = set(m["chain_id"] for m in receptor_metas)
    if len(chains) == 1:
        messages.append(f"PASS: All states use chain {chains.pop()}")
    else:
        messages.append(f"FAIL: Inconsistent chain IDs: {chains}")

    # Check residue range consistency
    ranges = [(m["state_name"], m["residue_start"], m["residue_end"], m["n_residues"])
              for m in receptor_metas]

    starts = set(r[1] for r in ranges)
    ends = set(r[2] for r in ranges)
    counts = set(r[3] for r in ranges)

    if len(starts) == 1 and len(ends) == 1:
        messages.append(f"PASS: All states share range {starts.pop()}-{ends.pop()}")
    else:
        messages.append(f"WARNING: Residue ranges differ:")
        for name, s, e, n in ranges:
            messages.append(f"  {name}: {s}-{e} ({n} residues)")

    if len(counts) == 1:
        messages.append(f"PASS: All states have {counts.pop()} residues")
    else:
        messages.append(f"WARNING: Residue counts differ:")
        for name, s, e, n in ranges:
            messages.append(f"  {name}: {n} residues")
        # Check if some states have gaps
        for m in receptor_metas:
            if m["missing_residues"]:
                messages.append(f"  {m['state_name']} missing: {m['missing_residues']}")

    # Check N-lobe presence
    for m in receptor_metas:
        if m["n_nlobe_residues"] > 0:
            messages.append(f"PASS: {m['state_name']} has N-lobe ({m['n_nlobe_residues']} res)")
        else:
            messages.append(f"FAIL: {m['state_name']} missing N-lobe!")

    # Check construct type
    types = set(m["construct_type"] for m in receptor_metas)
    if types == {"full_kinase_domain"}:
        messages.append(f"PASS: All states are full_kinase_domain construct")
    else:
        messages.append(f"WARNING: Mixed construct types: {types}")

    return messages


# ---------------------------------------------------------------------------
# Task 1.0.4: Partner validation
# ---------------------------------------------------------------------------

def validate_partner(partner_meta: dict) -> List[str]:
    """Validate extended beta-meander partner.

    Returns list of validation messages.
    """
    messages = []
    messages.append("")
    messages.append("=" * 60)
    messages.append("PARTNER VALIDATION")
    messages.append("=" * 60)

    # VAL962 artifact check
    if partner_meta["val962_is_first"]:
        messages.append("FAIL: VAL962 is still the first residue (artifact NOT eliminated)")
    else:
        messages.append(f"PASS: First residue is {partner_meta['first_residue']} (VAL962 artifact eliminated)")

    # Extension check
    if partner_meta["residue_start"] <= 955:
        messages.append(f"PASS: Construct starts at {partner_meta['residue_start']} (<= 955)")
    else:
        messages.append(f"WARNING: Construct starts at {partner_meta['residue_start']} (> 955, extension may be insufficient)")

    # Sheet completeness
    sheets = partner_meta["sheet_annotations"]
    all_complete = True
    for sname, sinfo in sheets.items():
        if sinfo["complete"]:
            messages.append(f"PASS: {sname} complete ({sinfo['present_count']}/{sinfo['total_count']} residues)")
        else:
            messages.append(f"WARNING: {sname} incomplete ({sinfo['present_count']}/{sinfo['total_count']} residues)")
            all_complete = False

    if all_complete:
        messages.append("PASS: All sheets (8-12) are complete")

    # Missing residues
    if partner_meta["missing_residues"]:
        messages.append(f"WARNING: Missing residues in construct: {partner_meta['missing_residues']}")
    else:
        messages.append("PASS: No gaps in construct")

    return messages


# ---------------------------------------------------------------------------
# Task 1.0.5: Metadata persistence
# ---------------------------------------------------------------------------

def write_receptor_metadata(receptor_metas: List[dict], output_dir: Path) -> Path:
    """Write receptor_metadata.csv."""
    output_path = output_dir / "receptor_metadata.csv"
    fieldnames = [
        "state_name", "construct_type", "source_pdb", "source_chain",
        "description", "output_pdb", "chain_id",
        "residue_start", "residue_end", "n_residues", "n_atoms",
        "n_nlobe_residues", "n_clobe_residues", "nlobe_clobe_boundary",
        "missing_residues", "numbering_system", "membrane_proximal_excluded",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(receptor_metas)
    return output_path


def write_partner_metadata(partner_meta: dict, output_dir: Path) -> Path:
    """Write partner_metadata.csv."""
    output_path = output_dir / "partner_metadata.csv"

    # Flatten sheet annotations for CSV
    flat = {k: v for k, v in partner_meta.items() if k != "sheet_annotations"}
    flat["active_face_residues"] = str(partner_meta.get("active_face_residues", []))
    flat["structural_support_residues"] = str(partner_meta.get("structural_support_residues", []))

    for sname, sinfo in partner_meta.get("sheet_annotations", {}).items():
        flat[f"{sname}_residues"] = str(sinfo["residues"])
        flat[f"{sname}_role"] = sinfo["role"]
        flat[f"{sname}_complete"] = sinfo["complete"]

    fieldnames = list(flat.keys())
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(flat)
    return output_path


def write_docking_pair_metadata(pair_metas: List[dict], output_dir: Path) -> Path:
    """Write docking_pair_metadata.csv."""
    output_path = output_dir / "docking_pair_metadata.csv"
    fieldnames = [
        "state_name", "output_pdb",
        "receptor_chain", "receptor_range", "receptor_n_residues",
        "partner_chain", "partner_range", "partner_n_residues",
        "excluded_residues_A", "n_excluded", "total_atoms",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(pair_metas)
    return output_path


# ---------------------------------------------------------------------------
# Task 1.0.6: Pilot data registration
# ---------------------------------------------------------------------------

def register_pilot_data(output_dir: Path) -> Path:
    """Register existing pilot data as historical reference."""
    output_path = output_dir / "pilot_data_reference.csv"

    pilot_entries = [
        {
            "label": "EGFR_dimer_beta_meander",
            "construct_type": "legacy_clobe_fragment",
            "receptor_pdb": "input/PPI/prepared/EGFR_dimer_beta_meander.pdb",
            "receptor_info": "input/PPI/prepared/EGFR_dimer_beta_meander_info.json",
            "receptor_description": "EGFR dimer (A+B merged) × C-lobe fragment (45 res)",
            "partner_pdb": "input/PPI/beta_meander.pdb",
            "partner_range": "960-1006",
            "known_artifacts": "N-lobe absence; VAL962 N-terminal artifact; no orientation filter",
            "status": "historical_reference_only",
            "results_dir": "output/egfr_myo1d_vina/ppi/beta_meander",
        },
        {
            "label": "EGFR_dimer_TH1",
            "construct_type": "legacy_clobe_fragment",
            "receptor_pdb": "input/PPI/prepared/EGFR_dimer_TH1.pdb",
            "receptor_info": "input/PPI/prepared/EGFR_dimer_TH1_info.json",
            "receptor_description": "EGFR dimer (A+B merged) × TH1 domain (206 res)",
            "partner_pdb": "input/PPI/TH1 domain.pdb",
            "partner_range": "801-1006",
            "known_artifacts": "N-lobe absence; TH1 noise (not primary search input)",
            "status": "historical_reference_only",
            "results_dir": "output/egfr_myo1d_vina/ppi/TH1",
        },
    ]

    if not pilot_entries:
        print("  WARNING: No pilot entries to register")
        return output_path
    fieldnames = list(pilot_entries[0].keys())
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(pilot_entries)

    print(f"\nPilot data registered: {output_path}")
    print(f"  {len(pilot_entries)} entries (historical reference only, not modified)")

    return output_path


# ---------------------------------------------------------------------------
# Validation report
# ---------------------------------------------------------------------------

def write_validation_report(
    receptor_metas: List[dict],
    partner_meta: dict,
    pair_metas: List[dict],
    receptor_messages: List[str],
    partner_messages: List[str],
    output_dir: Path,
) -> Path:
    """Write phase1_input_validation_report.md."""
    output_path = output_dir / "phase1_input_validation_report.md"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "# Phase 1 Input Validation Report",
        "",
        f"Generated: {timestamp}",
        f"Task Group: 1.0 (Receptor and Partner Input Preparation)",
        "",
        "## Summary",
        "",
        f"- **Receptor states prepared:** {len(receptor_metas)}",
        f"- **Receptor construct type:** full_kinase_domain",
        f"- **Partner construct type:** extended_beta_meander",
        f"- **Partner residue range:** {partner_meta['residue_start']}-{partner_meta['residue_end']}",
        f"- **Numbering system:** PDB (3GT8-consistent)",
        "",
        "## Receptor Details",
        "",
        "| State | Range | Residues | N-lobe | C-lobe | Gaps |",
        "|-------|-------|----------|--------|--------|------|",
    ]

    for m in receptor_metas:
        gaps = m["missing_residues"] if m["missing_residues"] else "None"
        lines.append(
            f"| {m['state_name']} | {m['residue_start']}-{m['residue_end']} | "
            f"{m['n_residues']} | {m['n_nlobe_residues']} | {m['n_clobe_residues']} | {gaps} |"
        )

    lines.extend([
        "",
        "## Partner Details",
        "",
        f"- **Construct:** Extended beta-meander",
        f"- **Source:** TH1 domain structure",
        f"- **Range:** {partner_meta['residue_start']}-{partner_meta['residue_end']} "
        f"({partner_meta['n_residues']} residues)",
        f"- **First residue:** {partner_meta['first_residue']}",
        f"- **VAL962 artifact eliminated:** {'Yes' if not partner_meta['val962_is_first'] else 'NO — REQUIRES ATTENTION'}",
        "",
        "### Sheet Annotations",
        "",
        "| Sheet | Residues | Role | Complete |",
        "|-------|----------|------|----------|",
    ])

    for sname, sinfo in partner_meta.get("sheet_annotations", {}).items():
        status = "Yes" if sinfo["complete"] else f"No ({sinfo['present_count']}/{sinfo['total_count']})"
        lines.append(f"| {sname} | {sinfo['residues']} | {sinfo['role']} | {status} |")

    lines.extend([
        "",
        "## Docking Pairs",
        "",
        "| State | Receptor | Partner | Excluded |",
        "|-------|----------|---------|----------|",
    ])

    for p in pair_metas:
        lines.append(
            f"| {p['state_name']} | A:{p['receptor_range']} ({p['receptor_n_residues']} res) | "
            f"B:{p['partner_range']} ({p['partner_n_residues']} res) | {p['n_excluded']} res |"
        )

    lines.extend([
        "",
        "## Validation Checks",
        "",
        "```",
    ])
    lines.extend(receptor_messages)
    lines.extend(partner_messages)
    lines.extend([
        "```",
        "",
        "## Numbering Reference",
        "",
        "- All residue numbers in this report use PDB numbering (3GT8-consistent)",
        "- UniProt (P00533) = PDB + 24 (approximate)",
        "- N-lobe/C-lobe boundary: residue 838 (PDB numbering)",
        "- Membrane-proximal excluded residues: see docking_pair_metadata.csv",
        "",
        "## Pilot Data",
        "",
        "- Legacy pilot data (C-lobe fragment + truncated beta-meander) registered as historical reference",
        "- See pilot_data_reference.csv for details",
        "- Pilot data is NOT used as validation target for new results",
        "",
    ])

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return output_path


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def prepare_phase1_inputs(output_dir: Optional[Path] = None) -> dict:
    """Prepare runtime Phase 1 docking inputs from source receptor/partner files."""
    output_dir = output_dir or PHASE1_RUNTIME_INPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Phase 1 — Task Group 1.0: Input Preparation")
    print(f"Output directory: {output_dir}")

    # -----------------------------------------------------------------------
    # Step 1: Prepare receptors
    # -----------------------------------------------------------------------
    receptor_metas = []
    for state_name, state_info in RECEPTOR_STATES.items():
        meta = prepare_receptor(state_name, state_info, output_dir)
        receptor_metas.append(meta)

    # -----------------------------------------------------------------------
    # Step 2: Prepare extended beta-meander
    # -----------------------------------------------------------------------
    partner_meta = prepare_extended_beta_meander(output_dir)

    # -----------------------------------------------------------------------
    # Step 3: Create docking pairs
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"Creating docking-ready PDB pairs")
    pair_metas = []
    for r_meta in receptor_metas:
        state_name = r_meta["state_name"]
        receptor_pdb = PROJECT_ROOT / r_meta["output_pdb"]
        partner_pdb = output_dir / "partner_extended_beta_meander.pdb"
        docking_pdb = output_dir / f"docking_{state_name}_ext_beta_meander.pdb"

        pair_meta = prepare_docking_pair(
            receptor_pdb, partner_pdb, docking_pdb, state_name
        )
        pair_metas.append(pair_meta)

    # -----------------------------------------------------------------------
    # Step 4: Validate
    # -----------------------------------------------------------------------
    receptor_messages = validate_receptors(receptor_metas)
    partner_messages = validate_partner(partner_meta)

    print("\n" + "\n".join(receptor_messages))
    print("\n".join(partner_messages))

    # -----------------------------------------------------------------------
    # Step 5: Write metadata
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("Writing metadata files")

    r_csv = write_receptor_metadata(receptor_metas, output_dir)
    print(f"  receptor_metadata.csv: {r_csv}")

    p_csv = write_partner_metadata(partner_meta, output_dir)
    print(f"  partner_metadata.csv: {p_csv}")

    d_csv = write_docking_pair_metadata(pair_metas, output_dir)
    print(f"  docking_pair_metadata.csv: {d_csv}")

    # -----------------------------------------------------------------------
    # Step 6: Register pilot data
    # -----------------------------------------------------------------------
    pilot_csv = register_pilot_data(output_dir)

    # -----------------------------------------------------------------------
    # Step 7: Write validation report
    # -----------------------------------------------------------------------
    report_path = write_validation_report(
        receptor_metas, partner_meta, pair_metas,
        receptor_messages, partner_messages, output_dir,
    )
    print(f"\nValidation report: {report_path}")

    # -----------------------------------------------------------------------
    # Final summary
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("TASK GROUP 1.0 COMPLETE")
    print(f"{'='*60}")
    print(f"\nDeliverables in {output_dir}/:")
    print(f"  Receptor PDBs:     3 states × full_kinase_domain")
    print(f"  Partner PDB:       1 × extended_beta_meander ({partner_meta['residue_start']}-{partner_meta['residue_end']})")
    print(f"  Docking pairs:     3 × (receptor + partner)")
    print(f"  receptor_metadata.csv")
    print(f"  partner_metadata.csv")
    print(f"  docking_pair_metadata.csv")
    print(f"  pilot_data_reference.csv")
    print(f"  phase1_input_validation_report.md")

    # Count pass/fail
    all_messages = receptor_messages + partner_messages
    n_pass = sum(1 for m in all_messages if m.startswith("PASS"))
    n_warn = sum(1 for m in all_messages if m.startswith("WARNING"))
    n_fail = sum(1 for m in all_messages if m.startswith("FAIL"))
    print(f"\nValidation: {n_pass} PASS, {n_warn} WARNING, {n_fail} FAIL")

    return {
        "output_dir": output_dir,
        "receptor_metadata_csv": r_csv,
        "partner_metadata_csv": p_csv,
        "docking_pair_metadata_csv": d_csv,
        "pilot_data_reference_csv": pilot_csv,
        "validation_report_path": report_path,
        "docking_pairs": {
            meta["state_name"]: output_dir / f"docking_{meta['state_name']}_ext_beta_meander.pdb"
            for meta in pair_metas
        },
        "status_code": 0 if n_fail == 0 else 1,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.0: Prepare full kinase domain + extended beta-meander inputs"
    )
    parser.add_argument(
        "--output_dir", type=Path,
        default=PHASE1_RUNTIME_INPUT_DIR,
        help="Output directory for prepared runtime inputs (default: output/phase1_ppi/runtime_inputs/)",
    )
    args = parser.parse_args()

    result = prepare_phase1_inputs(args.output_dir)
    return result["status_code"]


if __name__ == "__main__":
    sys.exit(main())
