#!/usr/bin/env python
"""Assess chemical diversity of the 3 project ligands.

Computes pairwise Morgan fingerprint Tanimoto similarity and
physicochemical descriptors for each ligand.  Outputs a CSV file
with pairwise comparisons and a diversity judgment.

Usage::

    python scripts/assess_ligand_diversity.py

Output: ``output/ligand_diversity_assessment.csv``
"""

import csv
import sys
from itertools import combinations
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    sys.exit("ERROR: RDKit is required.  Install with: pip install rdkit")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LIGAND_DIR = PROJECT_ROOT / "input" / "ligands"
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_CSV = OUTPUT_DIR / "ligand_diversity_assessment.csv"

LIGAND_FILES = [
    "173940_ligand.sdf",
    "97806_ligand.sdf",
    "VAX-C12_0_ligand.sdf",
]

# Tanimoto thresholds for diversity judgment
DIVERSE_THRESHOLD = 0.4
SIMILAR_THRESHOLD = 0.7


def _ligand_id(filename: str) -> str:
    return filename.replace("_ligand.sdf", "")


def load_molecule(sdf_path: Path) -> Chem.Mol:
    """Load first molecule from SDF file."""
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
    for mol in supplier:
        if mol is not None:
            return mol
    raise ValueError(f"Could not parse molecule from {sdf_path}")


def compute_descriptors(mol: Chem.Mol) -> dict:
    """Compute physicochemical descriptors."""
    return {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
        "RingCount": Descriptors.RingCount(mol),
        "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3),
    }


def compute_tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    """Compute Morgan fingerprint Tanimoto similarity."""
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_a = gen.GetFingerprint(mol_a)
    fp_b = gen.GetFingerprint(mol_b)
    return round(DataStructs.TanimotoSimilarity(fp_a, fp_b), 4)


def diversity_judgment(tanimoto: float) -> str:
    if tanimoto < DIVERSE_THRESHOLD:
        return "diverse"
    elif tanimoto <= SIMILAR_THRESHOLD:
        return "moderate"
    else:
        return "similar"


def main():
    # Load molecules
    molecules = {}
    for fname in LIGAND_FILES:
        path = LIGAND_DIR / fname
        if not path.exists():
            sys.exit(f"ERROR: Ligand file not found: {path}")
        mol = load_molecule(path)
        molecules[_ligand_id(fname)] = mol
        print(f"  Loaded {_ligand_id(fname)} from {fname}")

    # Compute per-ligand descriptors
    descriptors = {}
    for lig_id, mol in molecules.items():
        descriptors[lig_id] = compute_descriptors(mol)

    # Compute pairwise Tanimoto
    lig_ids = list(molecules.keys())
    pairwise = []
    for id_a, id_b in combinations(lig_ids, 2):
        tanimoto = compute_tanimoto(molecules[id_a], molecules[id_b])
        pairwise.append({
            "ligand_a": id_a,
            "ligand_b": id_b,
            "tanimoto_similarity": tanimoto,
            "diversity_judgment": diversity_judgment(tanimoto),
        })

    # Write CSV
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    header = [
        "row_type",
        "ligand_a", "ligand_b",
        "tanimoto_similarity", "diversity_judgment",
        "MW", "LogP", "TPSA", "HBD", "HBA",
        "RotatableBonds", "RingCount", "FractionCSP3",
    ]
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()

        # Descriptor rows (one per ligand)
        for lig_id in lig_ids:
            row = {"row_type": "descriptor", "ligand_a": lig_id}
            row.update(descriptors[lig_id])
            writer.writerow(row)

        # Pairwise rows
        for pair in pairwise:
            pair["row_type"] = "pairwise"
            writer.writerow(pair)

    print(f"\n  Output: {OUTPUT_CSV}")
    print(f"  Pairwise results:")
    for p in pairwise:
        print(f"    {p['ligand_a']} vs {p['ligand_b']}: "
              f"Tanimoto={p['tanimoto_similarity']}, "
              f"judgment={p['diversity_judgment']}")


if __name__ == "__main__":
    main()
