#!/usr/bin/env python3
"""Analyze 3GT8 EGFR dimer PDB to estimate membrane-proximal residues.
Uses only standard library (no numpy)."""

import math
from collections import defaultdict

PDB_FILE = "/home/user/codex_ligand/smoke_test/input/original/3gt8_dimer_anp.pdb"

def vec_sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def vec_add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def vec_scale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def vec_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vec_norm(a):
    return math.sqrt(vec_dot(a, a))

def vec_normalize(a):
    n = vec_norm(a)
    return (a[0]/n, a[1]/n, a[2]/n)

def vec_dist(a, b):
    return vec_norm(vec_sub(a, b))

def centroid(coords):
    n = len(coords)
    sx = sum(c[0] for c in coords)
    sy = sum(c[1] for c in coords)
    sz = sum(c[2] for c in coords)
    return (sx/n, sy/n, sz/n)

def parse_atoms(pdb_file):
    chains = defaultdict(list)
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            chain_resnum = line[21:26].strip()
            if chain_resnum[0] in "AB" and chain_resnum[1:].isdigit():
                chain = chain_resnum[0]
                resnum = int(chain_resnum[1:])
            else:
                chain = line[21].strip()
                resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            chains[chain].append((resnum, resname, atomname, x, y, z))
    return chains

def get_ca_atoms(atoms):
    ca = {}
    for resnum, resname, atomname, x, y, z in atoms:
        if atomname == "CA":
            ca[resnum] = ((x, y, z), resname)
    sorted_resnums = sorted(ca.keys())
    return ca, sorted_resnums

def find_ranges(residues, gap=1):
    if not residues:
        return []
    residues = sorted(residues)
    ranges = []
    start = residues[0]
    end = residues[0]
    for r in residues[1:]:
        if r <= end + gap + 1:
            end = r
        else:
            ranges.append((start, end))
            start = r
            end = r
    ranges.append((start, end))
    return ranges

def ranges_str(ranges):
    return ",".join(f"{s}-{e}" if s != e else str(s) for s, e in ranges)

def main():
    chains = parse_atoms(PDB_FILE)
    print(f"Chains found: {sorted(chains.keys())}")
    for ch in sorted(chains.keys()):
        print(f"  Chain {ch}: {len(chains[ch])} atoms")

    for chain_id in sorted(chains.keys()):
        ca, resnums = get_ca_atoms(chains[chain_id])
        n_res = len(resnums)
        print(f"\n{'='*70}")
        print(f"CHAIN {chain_id}: {n_res} residues, range {resnums[0]}-{resnums[-1]}")
        print(f"{'='*70}")

        n_term = 10
        c_term = 10
        nterm_resnums = resnums[:n_term]
        cterm_resnums = resnums[-c_term:]

        nterm_centroid = centroid([ca[r][0] for r in nterm_resnums])
        cterm_centroid = centroid([ca[r][0] for r in cterm_resnums])

        print(f"\nN-terminal residues (first {n_term}):")
        for r in nterm_resnums:
            coord = ca[r][0]
            resname = ca[r][1]
            print(f"  {chain_id}{r} {resname:3s}  ({coord[0]:7.2f}, {coord[1]:7.2f}, {coord[2]:7.2f})")
        print(f"  Centroid: ({nterm_centroid[0]:.2f}, {nterm_centroid[1]:.2f}, {nterm_centroid[2]:.2f})")

        print(f"\nC-terminal residues (last {c_term}):")
        for r in cterm_resnums:
            coord = ca[r][0]
            resname = ca[r][1]
            print(f"  {chain_id}{r} {resname:3s}  ({coord[0]:7.2f}, {coord[1]:7.2f}, {coord[2]:7.2f})")
        print(f"  Centroid: ({cterm_centroid[0]:.2f}, {cterm_centroid[1]:.2f}, {cterm_centroid[2]:.2f})")

        # Membrane direction: C-term -> N-term
        membrane_vec = vec_sub(nterm_centroid, cterm_centroid)
        membrane_vec_norm = vec_normalize(membrane_vec)
        print(f"\nMembrane direction vector (C-term -> N-term): "
              f"({membrane_vec_norm[0]:.3f}, {membrane_vec_norm[1]:.3f}, {membrane_vec_norm[2]:.3f})")
        print(f"  Vector length: {vec_norm(membrane_vec):.1f} A")

        # Project all CA atoms onto membrane direction
        projections = {}
        for r in resnums:
            coord = ca[r][0]
            proj = vec_dot(vec_sub(coord, cterm_centroid), membrane_vec_norm)
            projections[r] = proj

        proj_values = list(projections.values())
        proj_min = min(proj_values)
        proj_max = max(proj_values)
        proj_range = proj_max - proj_min

        threshold_20 = proj_max - 0.20 * proj_range
        membrane_proximal_proj = sorted(
            [r for r in resnums if projections[r] >= threshold_20],
            key=lambda r: -projections[r]
        )

        print(f"\n--- Projection Analysis ---")
        print(f"Projection range: {proj_min:.1f} to {proj_max:.1f} (span: {proj_range:.1f} A)")
        print(f"Top 20% threshold (membrane-proximal): projection >= {threshold_20:.1f}")
        print(f"\nMembrane-proximal residues (top 20%, {len(membrane_proximal_proj)} residues):")
        for r in membrane_proximal_proj:
            resname = ca[r][1]
            print(f"  {chain_id}{r} {resname:3s}  proj={projections[r]:.1f}")

        # Alternative: residues within 15A of N-terminal centroid
        nterm_close = sorted(
            [r for r in resnums if vec_dist(ca[r][0], nterm_centroid) <= 15.0],
            key=lambda r: vec_dist(ca[r][0], nterm_centroid)
        )

        print(f"\n--- Distance to N-terminal centroid (<= 15 A, {len(nterm_close)} residues) ---")
        for r in nterm_close:
            dist = vec_dist(ca[r][0], nterm_centroid)
            resname = ca[r][1]
            print(f"  {chain_id}{r} {resname:3s}  dist={dist:.1f} A")

        membrane_set = set(membrane_proximal_proj) | set(nterm_close)
        membrane_sorted = sorted(membrane_set)

        print(f"\n--- UNION of both methods ({len(membrane_sorted)} residues) ---")
        print(f"Residue list: {membrane_sorted}")
        ranges = find_ranges(membrane_sorted, gap=1)
        print(f"Contiguous ranges (allowing gap of 1):")
        for s, e in ranges:
            print(f"  {s}-{e} ({e-s+1} residues)")
        print(f"\nSuggested excluded_residues_{chain_id} (membrane only) = {ranges_str(ranges)}")

    # Dimer interface analysis
    print(f"\n{'='*70}")
    print("DIMER INTERFACE ANALYSIS")
    print(f"{'='*70}")

    ca_A, resnums_A = get_ca_atoms(chains["A"])
    ca_B, resnums_B = get_ca_atoms(chains["B"])

    interface_dist = 10.0
    interface_A = set()
    interface_B = set()

    for ra in resnums_A:
        for rb in resnums_B:
            d = vec_dist(ca_A[ra][0], ca_B[rb][0])
            if d < interface_dist:
                interface_A.add(ra)
                interface_B.add(rb)

    for ch, iface, ca_dict in [("A", sorted(interface_A), ca_A), ("B", sorted(interface_B), ca_B)]:
        print(f"\nChain {ch} dimer interface residues (CA within {interface_dist}A, {len(iface)} residues):")
        for r in iface:
            resname = ca_dict[r][1]
            print(f"  {ch}{r} {resname}")
        if iface:
            ranges = find_ranges(iface, gap=1)
            print(f"  Ranges: {ranges_str(ranges)}")

    # Combined
    print(f"\n{'='*70}")
    print("COMBINED EXCLUSION ZONES (membrane-proximal + dimer interface)")
    print(f"{'='*70}")

    for chain_id in ["A", "B"]:
        ca_dict, resnums = get_ca_atoms(chains[chain_id])
        nterm_resnums = resnums[:10]
        cterm_resnums = resnums[-10:]
        nterm_centroid = centroid([ca_dict[r][0] for r in nterm_resnums])
        cterm_centroid = centroid([ca_dict[r][0] for r in cterm_resnums])
        membrane_vec = vec_sub(nterm_centroid, cterm_centroid)
        membrane_vec_norm = vec_normalize(membrane_vec)

        projections = {}
        for r in resnums:
            proj = vec_dot(vec_sub(ca_dict[r][0], cterm_centroid), membrane_vec_norm)
            projections[r] = proj

        proj_values = list(projections.values())
        threshold_20 = max(proj_values) - 0.20 * (max(proj_values) - min(proj_values))

        membrane_set = set(r for r in resnums if projections[r] >= threshold_20)
        nterm_close = set(r for r in resnums if vec_dist(ca_dict[r][0], nterm_centroid) <= 15.0)

        iface = interface_A if chain_id == "A" else interface_B

        combined = sorted(membrane_set | nterm_close | iface)
        ranges = find_ranges(combined, gap=2)

        print(f"\nChain {chain_id} combined exclusion ({len(combined)} residues):")
        print(f"  Membrane-proximal: {len(membrane_set | nterm_close)} residues")
        print(f"  Dimer interface: {len(iface)} residues")
        print(f"  excluded_residues_{chain_id} = {ranges_str(ranges)}")

if __name__ == "__main__":
    main()
