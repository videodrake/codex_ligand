#!/usr/bin/env python3
"""Prepare EGFR dimer + MYO1D domain PDBs for PyRosetta PPI docking.

Workflow:
  1. Load EGFR dimer PDB (chain A + B)
  2. Remove non-protein chains (ANP ligands etc.)
  3. Renumber chain B residues with +1000 offset
  4. Merge into single chain A
  5. Add MYO1D domain as chain B
  6. Save combined PDB + mapping table for post-docking restoration

Usage:
  python -m egfr_pipeline.ppi.prepare_dimer_pdb \\
      --dimer smoke_test/input/original/3gt8_dimer_anp.pdb \\
      --partner input/PPI/beta_meander.pdb \\
      --output input/PPI/prepared/EGFR_dimer_beta_meander.pdb

Post-docking restoration:
  python -m egfr_pipeline.ppi.prepare_dimer_pdb --restore \\
      --input result.pdb \\
      --mapping input/PPI/prepared/EGFR_dimer_beta_meander_mapping.csv \\
      --output result_restored.pdb
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHAIN_B_OFFSET = 1000  # Chain B residues get +1000

# Membrane-proximal residues (user-provided from structure inspection)
# These residues face the membrane and should be excluded from PPI docking
MEMBRANE_PROXIMAL_A = "709-720,724-731,736-739,747,783-785,799-805,871-873,917-921"
MEMBRANE_PROXIMAL_B = "713-720,726-729,799-804,868-874,917-920"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_ranges(range_str: str) -> set:
    """Parse '699-711,748-770' into {699,700,...,711,748,...,770}."""
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


def format_ranges(nums: set) -> str:
    """Format {699,700,701,748,749} into '699-701,748-749'."""
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


def parse_pdb_atoms(pdb_path: Path) -> List[str]:
    """Read PDB and return ATOM/HETATM lines."""
    lines = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                lines.append(line)
    return lines


def get_chain(line: str) -> str:
    return line[21]


def get_resnum(line: str) -> int:
    return int(line[22:26].strip())


def set_chain(line: str, chain: str) -> str:
    return line[:21] + chain + line[22:]


def set_resnum(line: str, resnum: int) -> str:
    return line[:22] + f"{resnum:4d}" + line[26:]


def set_atom_serial(line: str, serial: int) -> str:
    return line[:6] + f"{serial:5d}" + line[11:]


def is_protein_atom(line: str) -> bool:
    """Check if line is a protein ATOM (not HETATM for ligands)."""
    return line.startswith("ATOM")


# ---------------------------------------------------------------------------
# Prepare: merge dimer + partner
# ---------------------------------------------------------------------------

def prepare_dimer_partner(
    dimer_path: Path,
    partner_path: Path,
    output_path: Path,
    partner_name: str = "",
) -> Tuple[Path, Path, dict]:
    """Merge EGFR dimer into chain A, add partner as chain B.

    Returns (output_pdb_path, mapping_csv_path, info_dict).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # --- Parse dimer ---
    dimer_lines = parse_pdb_atoms(dimer_path)
    chain_a_lines = [l for l in dimer_lines if is_protein_atom(l) and get_chain(l) == "A"]
    chain_b_lines = [l for l in dimer_lines if is_protein_atom(l) and get_chain(l) == "B"]

    # Residue ranges
    a_resnums = sorted(set(get_resnum(l) for l in chain_a_lines))
    b_resnums = sorted(set(get_resnum(l) for l in chain_b_lines))

    print(f"EGFR dimer: Chain A {a_resnums[0]}-{a_resnums[-1]} ({len(a_resnums)} res), "
          f"Chain B {b_resnums[0]}-{b_resnums[-1]} ({len(b_resnums)} res)")

    # --- Renumber chain B → chain A with offset ---
    merged_lines = list(chain_a_lines)  # chain A as-is
    mapping = {}  # new_resnum -> (original_chain, original_resnum)

    for r in a_resnums:
        mapping[r] = ("A", r)

    for line in chain_b_lines:
        old_resnum = get_resnum(line)
        new_resnum = old_resnum + CHAIN_B_OFFSET
        new_line = set_chain(line, "A")
        new_line = set_resnum(new_line, new_resnum)
        merged_lines.append(new_line)
        mapping[new_resnum] = ("B", old_resnum)

    new_b_resnums = sorted(set(get_resnum(set_resnum(l, get_resnum(l) + CHAIN_B_OFFSET))
                               for l in chain_b_lines))
    # Simpler: just compute
    new_b_range = (b_resnums[0] + CHAIN_B_OFFSET, b_resnums[-1] + CHAIN_B_OFFSET)

    print(f"Renumbered: Chain B {b_resnums[0]}-{b_resnums[-1]} → "
          f"A:{new_b_range[0]}-{new_b_range[1]} (offset +{CHAIN_B_OFFSET})")

    # --- Parse partner ---
    partner_lines = parse_pdb_atoms(partner_path)
    partner_protein = [l for l in partner_lines if is_protein_atom(l)]

    # Set partner to chain B
    partner_out = []
    for line in partner_protein:
        partner_out.append(set_chain(line, "B"))

    partner_resnums = sorted(set(get_resnum(l) for l in partner_protein))
    print(f"Partner ({partner_path.name}): {partner_resnums[0]}-{partner_resnums[-1]} "
          f"({len(partner_resnums)} res) → Chain B")

    # --- Combine and renumber atom serials ---
    all_lines = merged_lines + partner_out
    output_lines = []
    for i, line in enumerate(all_lines, 1):
        output_lines.append(set_atom_serial(line, i))

    # Write PDB
    with open(output_path, "w", encoding="utf-8") as f:
        for line in output_lines:
            f.write(line)
        f.write("TER\nEND\n")

    total_receptor = len(a_resnums) + len(b_resnums)
    print(f"Output: {output_path} ({len(output_lines)} atoms, "
          f"receptor={total_receptor} res, partner={len(partner_resnums)} res)")

    # --- Write mapping CSV ---
    mapping_path = output_path.with_name(output_path.stem + "_mapping.csv")
    with open(mapping_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["new_chain", "new_resnum", "original_chain", "original_resnum"])
        for new_r in sorted(mapping):
            orig_chain, orig_r = mapping[new_r]
            writer.writerow(["A", new_r, orig_chain, orig_r])
        for r in partner_resnums:
            writer.writerow(["B", r, "partner", r])

    print(f"Mapping: {mapping_path} ({len(mapping) + len(partner_resnums)} entries)")

    # --- Compute excluded residues ---
    excl_a_set = parse_ranges(MEMBRANE_PROXIMAL_A)
    excl_b_set = parse_ranges(MEMBRANE_PROXIMAL_B)
    # Offset chain B exclusion
    excl_b_offset = {r + CHAIN_B_OFFSET for r in excl_b_set}
    combined_excl = excl_a_set | excl_b_offset
    excl_str = format_ranges(combined_excl)

    info = {
        "dimer_path": str(dimer_path),
        "partner_path": str(partner_path),
        "partner_name": partner_name or partner_path.stem,
        "output_pdb": str(output_path),
        "mapping_csv": str(mapping_path),
        "chain_a_original": f"{a_resnums[0]}-{a_resnums[-1]}",
        "chain_b_original": f"{b_resnums[0]}-{b_resnums[-1]}",
        "chain_b_renumbered": f"{new_b_range[0]}-{new_b_range[1]}",
        "offset": CHAIN_B_OFFSET,
        "partner_residues": f"{partner_resnums[0]}-{partner_resnums[-1]}",
        "excluded_residues_A": excl_str,
        "n_excluded": len(combined_excl),
        "membrane_proximal_chainA": MEMBRANE_PROXIMAL_A,
        "membrane_proximal_chainB": MEMBRANE_PROXIMAL_B,
    }

    # Write info JSON
    info_path = output_path.with_name(output_path.stem + "_info.json")
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)

    print(f"\nExcluded residues (membrane-proximal, {len(combined_excl)} res):")
    print(f"  {excl_str}")
    print(f"\nInfo: {info_path}")

    return output_path, mapping_path, info


# ---------------------------------------------------------------------------
# Restore: post-docking chain renumbering
# ---------------------------------------------------------------------------

def restore_chains(
    input_pdb: Path,
    mapping_csv: Path,
    output_pdb: Path,
) -> Path:
    """Restore original chain IDs and residue numbers after docking."""
    # Load mapping
    chain_map: Dict[Tuple[str, int], Tuple[str, int]] = {}
    with open(mapping_csv, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            new_chain = row["new_chain"]
            new_resnum = int(row["new_resnum"])
            orig_chain = row["original_chain"]
            orig_resnum = int(row["original_resnum"])
            if orig_chain == "partner":
                # Partner stays as-is but gets a designated chain
                chain_map[(new_chain, new_resnum)] = ("C", orig_resnum)
            else:
                chain_map[(new_chain, new_resnum)] = (orig_chain, orig_resnum)

    # Process PDB
    lines = []
    with open(input_pdb, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = get_chain(line)
                resnum = get_resnum(line)
                key = (chain, resnum)
                if key in chain_map:
                    new_chain, new_resnum = chain_map[key]
                    line = set_chain(line, new_chain)
                    line = set_resnum(line, new_resnum)
                lines.append(line)
            elif line.startswith("TER") or line.startswith("END"):
                lines.append(line)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line)

    print(f"Restored: {output_pdb}")
    print(f"  Chain A = EGFR monomer 1 (original)")
    print(f"  Chain B = EGFR monomer 2 (original)")
    print(f"  Chain C = MYO1D partner")

    return output_pdb


def restore_csv(
    input_csv: Path,
    mapping_csv: Path,
    output_csv: Path,
    residue_columns: Optional[List[str]] = None,
) -> Path:
    """Restore residue numbering in a results CSV.

    Converts renumbered residue references back to original chain:resnum format.
    """
    # Load mapping: new_resnum -> (orig_chain, orig_resnum)
    resnum_map: Dict[int, Tuple[str, int]] = {}
    with open(mapping_csv, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["new_chain"] == "A":
                new_r = int(row["new_resnum"])
                orig_chain = row["original_chain"]
                orig_r = int(row["original_resnum"])
                if orig_chain != "partner":
                    resnum_map[new_r] = (orig_chain, orig_r)

    if residue_columns is None:
        residue_columns = ["binding_residues_A", "excluded_contacts_detail"]

    with open(input_csv, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    def restore_residue_ref(text: str) -> str:
        """Convert 'ALA1750' back to 'B:ALA750'."""
        def _replace(m):
            resname = m.group(1)
            resnum = int(m.group(2))
            if resnum in resnum_map:
                orig_chain, orig_r = resnum_map[resnum]
                return f"{orig_chain}:{resname}{orig_r}"
            return m.group(0)
        return re.sub(r'([A-Z]{3})(\d+)', _replace, text)

    for row in rows:
        for col in residue_columns:
            if col in row and row[col]:
                row[col] = restore_residue_ref(row[col])

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Restored CSV: {output_csv}")
    return output_csv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Prepare/restore EGFR dimer + MYO1D PDBs for PPI docking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command")

    # Prepare
    prep = sub.add_parser("prepare", help="Merge dimer + partner into docking-ready PDB")
    prep.add_argument("--dimer", required=True, help="EGFR dimer PDB path")
    prep.add_argument("--partner", required=True, help="MYO1D domain PDB path")
    prep.add_argument("--output", required=True, help="Output PDB path")
    prep.add_argument("--name", default="", help="Partner name (for labeling)")

    # Restore PDB
    rest = sub.add_parser("restore", help="Restore original chain numbering")
    rest.add_argument("--input", required=True, help="Docked result PDB")
    rest.add_argument("--mapping", required=True, help="Mapping CSV from prepare step")
    rest.add_argument("--output", required=True, help="Output restored PDB")

    # Restore CSV
    restcsv = sub.add_parser("restore-csv", help="Restore residue numbering in CSV")
    restcsv.add_argument("--input", required=True, help="Input CSV")
    restcsv.add_argument("--mapping", required=True, help="Mapping CSV from prepare step")
    restcsv.add_argument("--output", required=True, help="Output restored CSV")
    restcsv.add_argument("--columns", nargs="+", default=None,
                         help="Columns containing residue references")

    args = parser.parse_args()

    if args.command == "prepare":
        prepare_dimer_partner(
            Path(args.dimer), Path(args.partner), Path(args.output),
            partner_name=args.name,
        )
    elif args.command == "restore":
        restore_chains(
            Path(getattr(args, "input")), Path(args.mapping), Path(args.output),
        )
    elif args.command == "restore-csv":
        restore_csv(
            Path(getattr(args, "input")), Path(args.mapping), Path(args.output),
            residue_columns=args.columns,
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
